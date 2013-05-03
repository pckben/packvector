OBJ = main.o protein.o
OUT = out

GTEST = -L/opt/local/lib -lgtest -I/opt/local/include

all: test

$(OUT): $(OBJ)
	g++ -o $(OUT) $(OBJ) $(GTEST)

%.o: %.cc
	g++ -c $< -o $@ $(GTEST)

test: $(OUT)
	./$(OUT)

clean:
	-rm -f $(OBJ) $(OUT)
