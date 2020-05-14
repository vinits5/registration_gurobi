CCC = g++ -std=c++11
FLAGS = -g -Wall -Wextra 
GRBPATH = "/home/vinit/gurobi/gurobi752/linux64"
TREE = KDTree/kd_tree.h KDTree/kd_tree.cpp
DATA_TYPE = KDTree/type_defs.h
file = "main.cpp"
HELPER = helper.h helper.cpp
GUROBI_HELP = GurobiCodes/gurobi_helper.cpp GurobiCodes/gurobi_helper.h
ICP = ICP/ICP.cpp ICP/ICP.h
CALLBACK = Callback/callbackMtoS.h Callback/callbackMtoS.cpp
all:
	make compile && make run

compile:
	g++ -std=c++11 $(FLAGS) $(file) $(TREE) $(DATA_TYPE) $(HELPER) $(GUROBI_HELP) $(ICP) $(CALLBACK) -o exec -I $(GRBPATH)/include -L$(GRBPATH)/lib -lgurobi_c++ -lgurobi75

run:
	./exec

clean: 
	rm exec
