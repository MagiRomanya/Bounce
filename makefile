##
# Bounce
#
# @file
# @version 0.1

run:
	cc -o bounce main.c -lraylib -lm && ./bounce

test:
	cc -o test test.c && ./test


# end
