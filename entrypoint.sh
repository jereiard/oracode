#!/bin/bash
echo ">>> Starting Redis Service"
/usr/bin/redis-server --daemonize yes
echo ""
echo ">>> Starting Oracode"
echo "  Arguments: $@"
python /opt/reflower/reflow.py "$@"
