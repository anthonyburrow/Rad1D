CXX = g++
CXXFLAGS = -g

FC = gfortran
FFLAGS =

BIN_DIR = ./bin
SRC_DIR = ./src

# Directories
MAKE_DIRS = $(BIN_DIR) ./doc/figs

MKDIR_P = mkdir -p

.PHONY: directories $(MAKE_DIRS)

directories: $(MAKE_DIRS)

$(MAKE_DIRS):
	@$(MKDIR_P) $@

# Target 1
NAME1 = Rad1D
TARGET1 = $(BIN_DIR)/$(NAME1)
_SRC1 = $(NAME1)
_SRC2 = main RadModel initialize process io
SRC1 = $(_SRC1) $(_SRC2)

_OBJ_LIST = $(addsuffix .o, $(SRC1))
OBJ_LIST1 = $(addprefix $(BIN_DIR)/, $(_OBJ_LIST))

$(TARGET1): $(OBJ_LIST1)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Object recipes
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
#	$(DIR_GUARD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# $(BIN_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
# 	$(CC) $(CCFLAGS) -c $< -o $@

# $(BIN_DIR)/%.o: $(SRC_DIR)/%.f90
# 	$(FC) $(FFLAGS) -c $< -o $@

all: directories $(TARGET1)

TESTS = $(TARGET1)

test: directories $(TESTS)

clean:
	$(RM) -r $(BIN_DIR)