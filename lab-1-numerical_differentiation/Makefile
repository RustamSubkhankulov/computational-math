EXEC := comput

BUILD_DIR := ./build
RES_DIR   := ./res

SRC_DIRS := ./src
INC_DIRS := ./inc

#--------------------------------------

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_FLAGS := $(addprefix -I,$(INC_DIRS)) 
DEP_FLAGS := -MMD -MP

CXXFLAGS += --std=c++20 -O2

#--------------------------------------

.PHONY: build run clean plots

plots: $(RES_DIR)/comput_results.txt
	./scripts/plot.py

$(RES_DIR)/comput_res.txt: run

run: build $(RES_DIR)
	$(BUILD_DIR)/$(EXEC)

$(RES_DIR): 
	- mkdir -p $(RES_DIR) 

#--------------------------------------

build: $(BUILD_DIR)/$(EXEC)

$(BUILD_DIR)/$(EXEC): $(OBJS)
	$(CXX) $(^) -o $(@) $(LDFLAGS)

$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $(@))
	$(CXX) $(INC_FLAGS) $(DEP_FLAGS) $(CXXFLAGS) -c $(<) -o $(@)

#--------------------------------------

clean:
	- rm -r $(BUILD_DIR) 
	- rm -r $(RES_DIR)

#--------------------------------------

-include $(DEPS)