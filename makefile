##
# Bounce
#
# @file
# @version 0.1

all:
	cc -o bounce main.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11

test:
	cc -o test test.c && ./test


# end
