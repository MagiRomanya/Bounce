##
# Bounce
#
# @file
# @version 0.1

all:
	cc -o bounce main.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 && ./bounce

legacy:
	cc -o bounce main_backup.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 && ./bounce


# end
