#!/usr/bin/env bash

# Ivy must be installed in e.g. $HOME/.ant/lib

if [[ ! -f apache-ivy-2.5.0-bin.tar.gz ]]; then
  wget https://downloads.apache.org//ant/ivy/2.5.0/apache-ivy-2.5.0-bin.tar.gz
  echo "Fetched ivy"
fi

if [[ ! -f apache-ivy-2.5.0/ivy-2.5.0.jar ]]; then
  tar -zxvf apache-ivy-2.5.0-bin.tar.gz apache-ivy-2.5.0/ivy-2.5.0.jar
  echo "Unpacked apache-ivy-2.5.0/ivy-2.5.0.jar"
fi

if [[ ! -f ${HOME}/.ant/lib/ivy-2.5.0.jar ]]; then
  mkdir -p "${HOME}/.ant/lib/"
  cp apache-ivy-2.5.0/ivy-2.5.0.jar ${HOME}/.ant/lib/
  echo "Copied jar file to ${HOME}/.ant/lib/"
fi

echo "Ivy is installed"