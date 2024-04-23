
// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

#include <string>

class FeedbackToNetwork
{
public:

    //all static varialbes need to be initialized outside the class definition
    static std::string FeedbackSystemName;
    static std::string FeedbackvariableName;
    static Order Feedbackvariable_order;


    static void
    setup_feedback(EquationSystems* equation_systems);

    static void
    initialize_feedback(EquationSystems* equation_systems);

    static void
    assemble_feedback(EquationSystems& es,const std::string&  system_name);

    static void
    solve_feedback(EquationSystems* es);



  private:
      FeedbackToNetwork();
      FeedbackToNetwork(FeedbackToNetwork&);
      ~FeedbackToNetwork();
      FeedbackToNetwork& operator=(FeedbackToNetwork&);


  };
