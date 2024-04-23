// All rights reserved.
//

#ifndef included_ConoraryFlowInterface
#define included_ConoraryFlowInterface

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>


#include <string>

extern double   total_subdomain_volumes[100];
extern double   flow_change[100];
extern int  node_current_subdomain_lists[16403], element_current_subdomain_lists[16403], node_current_aha_lists[82290];

extern double total_flow_in, total_flow_out, total_flow_change, total_flow_difference, total_flow_difference_averaged_over_time;

extern double total_flow_difference_averaged_over_time_sub1;
extern double total_flow_difference_averaged_over_time_sub2;
extern double total_flow_difference_averaged_over_time_sub3;
extern double total_flow_difference_averaged_over_time_sub4;
extern double total_flow_difference_averaged_over_time_sub5;
extern double total_flow_difference_averaged_over_time_sub6;
extern double total_flow_difference_averaged_over_time_sub7;
extern double total_flow_difference_averaged_over_time_sub8;
extern double total_flow_difference_averaged_over_time_sub9;

class ConoraryFlowInterface
{
public:
    static bool send_source_from_coro;
    static bool visualise_coro_data;
    static bool output_coro_data_to_screen;
    static bool output_coro_data_to_file;
    static bool aha_to_subdomain;
    static bool exp_pao_to_file;

    static std::vector<int>  node_current_aha_lists;
    static std::vector<int>  node_current_subdomain_lists;
    static std::vector<int>  element_current_subdomain_lists;
    static double total_subdomain_volumes[100];
    static double flow_change[100];

    static double total_flow_in, total_flow_change, total_flow_out, total_flow_difference, total_flow_difference_averaged_over_time;

    static double total_flow_difference_averaged_over_time_sub1;
    static double total_flow_difference_averaged_over_time_sub2;
    static double total_flow_difference_averaged_over_time_sub3;
    static double total_flow_difference_averaged_over_time_sub4;
    static double total_flow_difference_averaged_over_time_sub5;
    static double total_flow_difference_averaged_over_time_sub6;
    static double total_flow_difference_averaged_over_time_sub7;
    static double total_flow_difference_averaged_over_time_sub8;
    static double total_flow_difference_averaged_over_time_sub9;


    //all static varialbes need to be initialized outside the class definition
    static std::string Pext_systemName;
    static std::string Pext_variableName;
    static Order Pext_variable_order;

    static std::string subdomain_test_systemName;
    static std::string subdomain_test_variableName;
    static Order subdomain_test_variable_order;

    static std::string visualise_coro_data_systemName;
    static std::string visualise_coro_data_variableName;
    static Order visualise_coro_data_variable_order;

    static std::string visualise_source_data_systemName;
    static std::string visualise_source_data_variableName;
    static Order visualise_source_data_variable_order;

    static std::string visualise_sink_data_systemName;
    static std::string visualise_sink_data_variableName;
    static Order visualise_sink_data_variable_order;


	static std::string Q_elem_systemName;
    static std::string Q_elem_variableName;
    static Order Q_elem_variable_order;


	static std::vector<double> totalVolume;

//---------------------------------------------------------
// setting up the P_ext system to interact with the coronary flow model
    static void
    setup_Pext(
        EquationSystems* equation_systems);

    static void
    initialize_Pext(
            EquationSystems* equation_systems);

   static void
     update_Pext(EquationSystems* equation_systems,double time);


    static void
    setup_Q_elem(
            EquationSystems* equation_systems);
    static void
    initialize_Q_elem(
            EquationSystems* equation_systems);

    static void
     update_Q_elem(EquationSystems* equation_systems,const int);

     static void
      update_Q_elem_hex(EquationSystems* equation_systems,const int);


    static void
    volume_of_subdomain(EquationSystems* equation_systems,
           std::vector<double>& volume_subdomain, unsigned int total_subdomain);



           static void
           setup_id_test(
               EquationSystems* equation_systems);

           static void
           initialize_id_test(
                   EquationSystems* equation_systems);

          static void
            update_id_test(EquationSystems* equation_systems);



            static void
            setup_visualise_coro_data(
                EquationSystems* equation_systems);

            static void
            initialize_visualise_coro_data(
                    EquationSystems* equation_systems);

           static void
             update_visualise_coro_data(EquationSystems* equation_systems);


             static void
               aha_to_subdomains(EquationSystems* equation_systems);

             static void
               return_Pext_to_zero(EquationSystems* equation_systems,double time,double dt);



//---------------------------------------------------------


private:
    ConoraryFlowInterface();
    ConoraryFlowInterface(ConoraryFlowInterface&);
    ~ConoraryFlowInterface();
    ConoraryFlowInterface& operator=(ConoraryFlowInterface&);
};


#endif
