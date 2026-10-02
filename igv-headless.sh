#!/bin/sh
set -eu

state_dir="$(mktemp -d)"
display_file="$state_dir/display"
xvfb_log="$state_dir/xvfb.log"
xvfb_pid=""

cleanup() {
    if [ -n "$xvfb_pid" ]; then
        kill "$xvfb_pid" 2>/dev/null || true
        wait "$xvfb_pid" 2>/dev/null || true
    fi
    rm -rf "$state_dir"
}
trap cleanup EXIT INT TERM

Xvfb -displayfd 3 -screen 0 1280x1024x24 -nolisten tcp \
    3>"$display_file" 2>"$xvfb_log" &
xvfb_pid=$!

attempt=0
while [ ! -s "$display_file" ]; do
    if ! kill -0 "$xvfb_pid" 2>/dev/null; then
        cat "$xvfb_log" >&2
        exit 1
    fi
    attempt=$((attempt + 1))
    if [ "$attempt" -ge 100 ]; then
        echo "Xvfb did not publish a display within 10 seconds" >&2
        cat "$xvfb_log" >&2
        exit 1
    fi
    sleep 0.1
done

DISPLAY=":$(cat "$display_file")"
export DISPLAY
"$@"
