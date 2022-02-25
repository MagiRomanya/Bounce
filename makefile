##
# Bounce
#
# @file
# @version 0.1

all:
	cc -o bounce main.c -lraylib -lm

test:
	cc -o test test.c && ./test


# end
