NAME = accretion_disk_simulator

CC = gcc

CFLAGS = -std=c99 -O3 -g -g3 -lm -O3 -lgmp -lGLU -lGL -lglut -fopenmp-simd -ffast-math -march=native -fpack-struct -fno-strict-aliasing -fno-omit-frame-pointer -fno-stack-protector -fno-exceptions 

SRC = $(wildcard *.c)

OBJ = $(SRC:.c=.o)

all: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(OBJ) -o $(NAME) $(CFLAGS)

clean:
	rm -f $(OBJ)

fclean: clean
	rm -f $(NAME)

re: fclean all

.PHONY: all clean fclean re