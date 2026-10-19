"""Provider adapters: run one question under one condition against one model.

Uniform contract -- every adapter returns:
    {
      "final_text": str,          # model's final assistant text
      "thinking":  str,           # extended-thinking / reasoning text, if exposed
      "tool_calls": [ {name, args, result} ],
      "raw": {...},               # finish reason, usage, iterations
      "error": str | None,
    }

Conditions:
    "unaugmented" -> no tools attached
    "surfmcp"     -> the 25 SurfactantKit tools attached via the provider's own
                     function-calling interface, executed against the real library

Only stdlib urllib is used, so the harness has no SDK version dependencies.
Keys come from environment variables only -- never a file in this (public) repo.
"""

from __future__ import annotations
import json
import os
import time
import urllib.request
import urllib.error

from tools import dispatch, to_openai_tools, to_anthropic_tools, to_gemini_tools

MAX_TOOL_ITERS = 8
MAX_RETRIES = 5


def _post(url, headers, body, timeout=300):
    data = json.dumps(body).encode()
    for attempt in range(MAX_RETRIES):
        req = urllib.request.Request(url, data=data, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 529) and attempt < MAX_RETRIES - 1:
                retry_after = e.headers.get("Retry-After")
                wait = float(retry_after) if retry_after else (2 ** attempt) * 2
                time.sleep(min(wait, 60))
                continue
            raise


# =============================== Anthropic-format (Claude + Qwen/DashScope) ===
def run_anthropic_like(model, question_system, question_user, condition,
                       base_url, api_key, thinking_budget=0):
    tools = to_anthropic_tools() if condition == "surfmcp" else None
    messages = [{"role": "user", "content": question_user}]
    tool_calls = []
    headers = {"x-api-key": api_key, "anthropic-version": "2023-06-01",
               "content-type": "application/json"}
    # max_tokens must cover thinking + the final response; thinking tokens count
    # against the same budget, so a fixed max_tokens smaller than thinking_budget
    # (or even just tight against a verbose unaugmented derivation) silently
    # truncates the answer before FINAL_JSON is ever reached.
    max_tokens = max(8192, thinking_budget + 4096)
    body_base = {"model": model, "max_tokens": max_tokens, "system": question_system}
    if thinking_budget:
        if "api.anthropic.com" in base_url:
            # Real Anthropic API changed its extended-thinking schema since this
            # harness was built: "thinking.type.enabled"/"budget_tokens" now
            # returns a 400 invalid_request_error on claude-opus-4-8 -- confirmed
            # directly against the live API (not guessed), error message names the
            # real replacement. DashScope's Anthropic-compatible endpoint (Qwen)
            # still accepts the old budget_tokens form (confirmed working in the
            # same pilot run), so this branch is Anthropic-API-specific, not
            # applied to the shared run_anthropic_like() path unconditionally.
            body_base["thinking"] = {"type": "adaptive"}
            body_base["output_config"] = {"effort": "high"}
        else:
            body_base["thinking"] = {"type": "enabled", "budget_tokens": thinking_budget}
    if tools:
        body_base["tools"] = tools

    thinking_text = []
    for _ in range(MAX_TOOL_ITERS):
        body = dict(body_base, messages=messages)
        resp = _post(base_url.rstrip("/") + "/v1/messages", headers, body)
        blocks = resp.get("content", [])
        for b in blocks:
            if b.get("type") == "thinking":
                thinking_text.append(b.get("thinking", ""))
        stop = resp.get("stop_reason")
        if stop == "tool_use":
            messages.append({"role": "assistant", "content": blocks})
            results = []
            for b in blocks:
                if b.get("type") == "tool_use":
                    out = dispatch(b["name"], b.get("input", {}))
                    tool_calls.append({"name": b["name"], "args": b.get("input", {}), "result": out})
                    results.append({"type": "tool_result", "tool_use_id": b["id"],
                                    "content": json.dumps(out)})
            messages.append({"role": "user", "content": results})
            continue
        final = "".join(b.get("text", "") for b in blocks if b.get("type") == "text")
        return {"final_text": final, "thinking": "\n".join(thinking_text),
                "tool_calls": tool_calls,
                "raw": {"stop_reason": stop, "usage": resp.get("usage")}, "error": None}
    return {"final_text": "", "thinking": "\n".join(thinking_text), "tool_calls": tool_calls,
            "raw": {"stop_reason": "max_iters"}, "error": "exceeded MAX_TOOL_ITERS"}


# =============================== OpenAI Chat Completions ======================
def run_openai(model, question_system, question_user, condition, api_key,
               reasoning_effort=None):
    tools = to_openai_tools() if condition == "surfmcp" else None
    messages = [{"role": "system", "content": question_system},
                {"role": "user", "content": question_user}]
    tool_calls = []
    headers = {"Authorization": f"Bearer {api_key}", "content-type": "application/json"}

    for _ in range(MAX_TOOL_ITERS):
        body = {"model": model, "messages": messages}
        if tools:
            body["tools"] = tools
            body["tool_choice"] = "auto"
        if reasoning_effort:
            body["reasoning_effort"] = reasoning_effort
        resp = _post("https://api.openai.com/v1/chat/completions", headers, body)
        msg = resp["choices"][0]["message"]
        if msg.get("tool_calls"):
            messages.append(msg)
            for tc in msg["tool_calls"]:
                name = tc["function"]["name"]
                try:
                    args = json.loads(tc["function"].get("arguments") or "{}")
                except Exception:
                    args = {}
                out = dispatch(name, args)
                tool_calls.append({"name": name, "args": args, "result": out})
                messages.append({"role": "tool", "tool_call_id": tc["id"],
                                 "content": json.dumps(out)})
            continue
        return {"final_text": msg.get("content") or "", "thinking": "",
                "tool_calls": tool_calls,
                "raw": {"finish_reason": resp["choices"][0].get("finish_reason"),
                        "usage": resp.get("usage")}, "error": None}
    return {"final_text": "", "thinking": "", "tool_calls": tool_calls,
            "raw": {"finish_reason": "max_iters"}, "error": "exceeded MAX_TOOL_ITERS"}


# =============================== Gemini generateContent =======================
def run_gemini(model, question_system, question_user, condition, api_key,
               thinking_budget=None, image_base64=None, image_mime_type="image/png"):
    tools = to_gemini_tools() if condition == "surfmcp" else None
    user_parts = [{"text": question_user}]
    if image_base64:
        # Real vision questions (GZ-15/36-40) require an attached chart -- Gemini's
        # REST API takes inline images as a base64 part alongside the text part,
        # both in the same user turn's "parts" list.
        user_parts.append({"inline_data": {"mime_type": image_mime_type, "data": image_base64}})
    contents = [{"role": "user", "parts": user_parts}]
    tool_calls = []
    url = (f"https://generativelanguage.googleapis.com/v1beta/models/{model}"
           f":generateContent?key={api_key}")
    headers = {"content-type": "application/json"}

    for _ in range(MAX_TOOL_ITERS):
        body = {"contents": contents,
                "systemInstruction": {"parts": [{"text": question_system}]}}
        if tools:
            body["tools"] = tools
        if thinking_budget is not None:
            body["generationConfig"] = {"thinkingConfig": {"thinkingBudget": thinking_budget}}
        resp = _post(url, headers, body)
        cand = resp["candidates"][0]
        parts = cand.get("content", {}).get("parts", [])
        fcalls = [p["functionCall"] for p in parts if "functionCall" in p]
        if fcalls:
            contents.append({"role": "model", "parts": parts})
            resp_parts = []
            for fc in fcalls:
                out = dispatch(fc["name"], dict(fc.get("args", {})))
                tool_calls.append({"name": fc["name"], "args": dict(fc.get("args", {})), "result": out})
                resp_parts.append({"functionResponse": {"name": fc["name"], "response": out}})
            contents.append({"role": "user", "parts": resp_parts})
            continue
        final = "".join(p.get("text", "") for p in parts if "text" in p)
        return {"final_text": final, "thinking": "", "tool_calls": tool_calls,
                "raw": {"finish_reason": cand.get("finishReason"),
                        "usage": resp.get("usageMetadata")}, "error": None}
    return {"final_text": "", "thinking": "", "tool_calls": tool_calls,
            "raw": {"finish_reason": "max_iters"}, "error": "exceeded MAX_TOOL_ITERS"}


# =============================== top-level dispatch ===========================
# model registry: logical name -> (provider, model_id, thinking-config-note)
def run(provider, model_id, system, user, condition, **kw):
    if provider == "anthropic":
        key = os.environ.get("ANTHROPIC_API_KEY")
        if not key:
            return _keyerr("ANTHROPIC_API_KEY")
        return run_anthropic_like(model_id, system, user, condition,
                                  "https://api.anthropic.com", key,
                                  thinking_budget=kw.get("thinking_budget", 0))
    if provider == "qwen":
        key = os.environ.get("DASHSCOPE_API_KEY")
        if not key:
            return _keyerr("DASHSCOPE_API_KEY")
        return run_anthropic_like(model_id, system, user, condition,
                                  "https://dashscope-intl.aliyuncs.com/apps/anthropic", key,
                                  thinking_budget=kw.get("thinking_budget", 0))
    if provider == "openai":
        key = os.environ.get("OPENAI_API_KEY")
        if not key:
            return _keyerr("OPENAI_API_KEY")
        return run_openai(model_id, system, user, condition, key,
                          reasoning_effort=kw.get("reasoning_effort"))
    if provider == "gemini":
        key = os.environ.get("GEMINI_API_KEY")
        if not key:
            return _keyerr("GEMINI_API_KEY")
        return run_gemini(model_id, system, user, condition, key,
                          thinking_budget=kw.get("thinking_budget"),
                          image_base64=kw.get("image_base64"),
                          image_mime_type=kw.get("image_mime_type", "image/png"))
    return {"final_text": "", "thinking": "", "tool_calls": [], "raw": {},
            "error": f"unknown provider {provider}"}


def _keyerr(var):
    return {"final_text": "", "thinking": "", "tool_calls": [], "raw": {},
            "error": f"missing env var {var}"}
