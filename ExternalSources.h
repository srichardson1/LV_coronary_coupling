// All rights reserved.
//

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

// ExternalSources is a static class that provides data and functions
// required to implement loading conditions and boundary constraints for the LV
// model.
class ExternalSources
{
   public:
    

    static void
    external_source_function(
        double& S,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const vector<const vector<double>*>& system_data,
        const vector<const vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
        double time,
        void* ctx);

  private:
    ExternalSources();
    ExternalSources(ExternalSources&);
    ~ExternalSources();
    ExternalSources& operator=(ExternalSources&);
};
