// APPLICATION INCLUDES
#include <ExternalSources.h>


// STATIC VARIABLES


// CLASS IMPLEMENTATION
void
ExternalSources::external_source_function(
                                          double& S,
                                          const TensorValue<double>& FF,
                                          const libMesh::Point& x,
                                          const libMesh::Point& X,
                                          Elem* const elem,
                                          const vector<const vector<double>*>& var_data,
                                          const vector<const vector<libMesh::VectorValue<double> >*>& grad_var_data,
                                          double time,
                                          void* ctx)
{


    const std::vector<double> & s_star = *var_data[0];

    //convert into reference state
    //const std::vector<double>& J_vec = *var_data[0];
    //const std::vector<libMesh::VectorValue<double> >& grad_J = *grad_var_data[0];

    double source = s_star[0];

    S = source;

//S = 1.0/65536;

    return;
}
