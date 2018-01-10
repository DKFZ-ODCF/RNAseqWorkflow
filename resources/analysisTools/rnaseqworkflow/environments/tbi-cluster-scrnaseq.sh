#!/usr/bin/env bash

source "${TOOL_WORKFLOW_ENVIRONMENT_TBI_CLUSTER}"
export UMITOOLS_BINARY=umi_tools

# module load "Je/${JE_VERSION:?No JE_VERSION}" # We won't use Je 1.0 since it appends extra space in read name
module load "Jemultiplexer/${JEMULTIPLEXER_VERSION:?No JEMULTIPLEXER_VERSION}"
export JE_BINARY=jemultiplexer