#!/bin/bash

vagrant up

vagrant ssh -c "cd /home/ubuntu/renjin && mvn -f vagrantpom.xml clean install"
