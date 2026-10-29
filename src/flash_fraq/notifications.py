"""Discord webhook notifications for pipeline job completion/failure.

Adapted from the `notify_bot` pattern in star-platinum/helpers.py, using
`requests` (already a project dependency) instead of `discord.py` since this
project has no other use for the discord client library.
"""

from __future__ import annotations

import os

import requests

# Read from the environment so no credential lives in the repo. Unset = notifications are a no-op.
DISCORD_WEBHOOK_ENV_VAR = "DISCORD_WEBHOOK_URL"


def notify_discord(msg: str) -> None:
    webhook_url = os.environ.get(DISCORD_WEBHOOK_ENV_VAR)
    if not webhook_url:
        return
    if len(msg) > 2000:
        msg = f"{msg[:1993]} [...]"
    requests.post(webhook_url, json={"content": msg}, timeout=10)
