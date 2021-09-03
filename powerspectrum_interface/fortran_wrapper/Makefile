ifeq (${FC},)
FC = gfortran
endif
ifeq (${FC},f77)
FC = gfortran
endif

# Source-code directory
LIB_SRC_DIR = ../../library/src

# Build directory
BUILD_DIR = build

LIB_DIR = lib
INCLUDE_DIR = include

# Debug build directory
DEBUG_BUILD_DIR = debug_build

ifeq (,$(findstring gfortran,${FC}))
FFLAGS = -fpic -fpp -qopenmp -fp-model precise -W0 -WB -O3 -ipo -axCORE-AVX2 -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow
DEBUGFLAGS = -fpp -qopenmp -g -check all -check noarg_temp_created -traceback -fpe0
MODOUT = -module $(INCLUDE_DIR)
DEBUG_MODOUT = -module $(INCLUDE_DIR)
else
FFLAGS =  -fPIC -cpp -fopenmp -O3 -ffree-line-length-none -std=gnu -fdefault-real-8 -fdefault-double-8
DEBUGFLAGS = -g -O0 -fcheck=all -fbacktrace -cpp -fdec -fopenmp -ffree-line-length-none -Wall -fbounds-check -ffpe-trap=invalid,zero,overflow
#add -fpe0 to check for floating point errors (think lowLike also throws these harmlessly)
MODOUT = -J$(INCLUDE_DIR)
DEBUG_MODOUT = -J$(INCLUDE_DIR)
endif

# Objects
# _OBJ = \
# 	precision.o \
# 	constants.o \
# 	physics.o \
# 	sorting.o \
# 	basic_operations.o \
# 	special_functions.o \
# 	array_operations.o \
# 	file_info.o \
# 	random_numbers.o \
# 	table_integer.o \
# 	interpolate.o \
# 	solve_equations.o \
# 	string_operations.o \
# 	calculus_table.o \
# 	statistics.o \
# 	calculus.o \
# 	minimization.o \
# 	camb_stuff.o \
# 	multidark_stuff.o \
# 	cosmology_functions.o \
# 	hmx.o \

_OBJ = \
	precision.o \
	constants.o \
	physics.o \
	sorting.o \
	basic_operations.o \
	special_functions.o \
	array_operations.o \
	file_info.o \
	random_numbers.o \
	table_integer.o \
	interpolate.o \
	ode_solvers.o \
	string_operations.o \
	calculus_table.o \
	statistics.o \
	calculus.o \
	minimization.o \
	root_finding.o \
	camb_stuff.o \
	multidark_stuff.o \
	hod_functions.o \
	cosmology_functions.o \
	hmx.o \


# Add prefixes of build directory to objects
OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))
DEBUG_OBJ = $(addprefix $(DEBUG_BUILD_DIR)/,$(_OBJ))

# ?
make_dirs = @mkdir -p $(@D)

all: wrapper

# Standard rules
lib: $(LIB_DIR)/libhmx.a
wrapper: $(LIB_DIR)/libhmx_wrapper.so

# Debugging rules
debug: FFLAGS += $(DEBUGFLAGS)
debug: $(LIB_DIR)/libhmx_wrapper.so

# Rule to make object files
$(BUILD_DIR)/%.o: $(LIB_SRC_DIR)/%.f90 $(INCLUDE_DIR)
	$(make_dirs)
	$(FC) -c -o $@ $< $(MODOUT) $(LDFLAGS) $(FFLAGS)

# Rule to make debugging objects
$(DEBUG_BUILD_DIR)/%.o: $(LIB_SRC_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< $(DEBUG_MODOUT) $(LDFLAGS) $(FFLAGS)

# Create include directory for .mod files
$(INCLUDE_DIR):
	@mkdir -p $(INCLUDE_DIR)
	
# Rule to make HMx static library
$(LIB_DIR)/libhmx.a: $(OBJ)
	@echo "\nBuilding static library.\n"
	$(make_dirs)
	$(AR) rc $@ $^

$(LIB_DIR)/libhmx_wrapper.so: HMx_wrapper.f90 $(LIB_DIR)/libhmx.a
	@echo "\nBuilding HMx wrapper.\n"
	$(make_dirs)
	$(FC) -o $@ $^ -I$(INCLUDE_DIR) $(LDFLAGS) -shared $(FFLAGS)

# Clean up
.PHONY: clean
clean:
	rm -f $(LIB_DIR)/libhmx.a
	rm -f $(LIB_DIR)/libhmx_wrapper.so
	test -n "$(LIB_DIR)" && rm -rf $(LIB_DIR)/libhmx_wrapper.so.dSYM/
	rm -f $(BUILD_DIR)/*.o
	rm -f $(INCLUDE_DIR)/*.mod
	rm -f $(SRC_DIR)/*.mod
	rm -f *.mod
	rm -f $(DEBUG_BUILD_DIR)/*.o

