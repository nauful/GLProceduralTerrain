APP_NAME = terrain

CC = g++ -DLINUX
LIBS = -lGL -lX11 -lglut -lGLEW
CC_OBJ = $(CC) -c

RELEASE = -O3 -Wno-unused-result -ffast-math

SRCS = $(GLEW_SRC) $(shell ls src/*.cpp) $(shell ls lib/*.cpp)

rel: $(SRCS)
	$(CC) $(RELEASE) $(SRCS) $(LIBS) -o build/$(APP_NAME)

all: rel

clean:
	rm build/$(APP_NAME)
