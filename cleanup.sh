#!/bin/bash

# 当前0级数据处理完成，告知下阶段继续处理
# touch /next-stage/$1
send-message $1
