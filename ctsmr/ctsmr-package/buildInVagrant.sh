#!/bin/bash

RED='\033[0;31m'
NC='\033[0m' # No Color

if ! command -v vagrant; then
  echo "${RED}vagrant is not available in path, is it installed?${NC}"
  exit
fi

vagrant up

vagrant ssh -c "cd /home/ubuntu/renjin && mvn -f vagrantpom.xml clean install"
