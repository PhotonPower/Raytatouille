# Runs EXE with the arguments ARGS (a ;-list) and fails unless its exit code is EXPECT.
# Used by the CLI tests that must tell exit 1 (input rejected) from exit 2 (wrong call).
execute_process(COMMAND ${EXE} ${ARGS} RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
if(NOT result EQUAL EXPECT)
  message(FATAL_ERROR "exit code ${result}, expected ${EXPECT}")
endif()
