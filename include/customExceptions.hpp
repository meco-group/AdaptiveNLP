#ifndef __CUSTOM_EXCEPTIONS__
#define __CUSTOM_EXCEPTIONS__
#include <exception>

// Exception thrown when there is an attempt to define building blocks when
// they are already defined before.
class illegalBlocksChangeException : public std::exception {
    public:
        const char * what () const throw() {
        return msg;
    }
    private:
        const char* msg = "Illegal attempt to change building blocks. Once "
                          "initialized, building blocks cannot change.";
};

class illegalDiscBlocksProvided : public std::exception {
    public:
        const char * what () const throw() {
        return msg;
    }
    private:
        const char* msg = "Illegal discretization function blocks provided. "
                          "All discretization blocks should have a different "
                          "number of time-steps as inputs.";
};

// Exception thrown when there is an attempt to remove an extra constraint
// that does not exist
class illegalConstraintRemovalException : public std::exception {
    public:
        const char * what () const throw() {
        return msg;
    }
    private:
        const char* msg = "Cannot remove unexisting extra constraint.";
};

// Expection thrown when a change is made to the NLP which is invalid
class illegalChangeToNLP : public std::exception {
    public:
        illegalChangeToNLP(std::string cause_msg){
            error_msg = 
                "Illegal change to the NLP problem structure (" + cause_msg + ")";
        }

        const char * what() const throw() {
            return error_msg.c_str();
        }
    private:
        std::string error_msg;
};

// Exception thrown when invalid parameters are provided
class illegalParameterProvided : public std::exception {
    public:

        const char * what() const throw() {
            return msg;
        }
    private:
        const char* msg = "Parameters provided do not have the right dimensions";
};

class unexistingParameterProvided : public std::exception {
    public:
        const char * what() const throw() {
            return msg;
        }
    private:
        const char* msg = "Discretization parameter provided for "
                          "discretization function that does not exist";
};

class illegalCallToSolve : public std::exception {
    public:
        const char * what() const throw() {
            return msg;
        }
    private:
        const char* msg = "A call to solveNLP is not allowed if the problem "
                          "structure has been changed. Please call "
                          "clearStructuralZeros() first";
};

#endif