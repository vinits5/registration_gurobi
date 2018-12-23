CCC = g++ -std=c++11
FLAGS = -g -Wall -Wextra 
GRBPATH = "/home/vinit/gurobi/gurobi752/linux64"
TREE = KDTree/kd_tree.h KDTree/kd_tree.cpp
DATA_TYPE = KDTree/type_defs.h
file = "tejas2.cpp"
HELPER = helper.h helper.cpp
GUROBI_HELP = GurobiCodes/gurobi_helper.cpp GurobiCodes/gurobi_helper.h
ICP = ICP/ICP.cpp ICP/ICP.h
all:
	make compile && make run

compile:
	g++ -std=c++11 $(FLAGS) $(file) $(TREE) $(DATA_TYPE) $(HELPER) $(GUROBI_HELP) $(ICP) -o exec -I $(GRBPATH)/include -L$(GRBPATH)/lib -lgurobi_c++ -lgurobi75

run:
	./exec

clean: 
	rm exec
