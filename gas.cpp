#include "atomic.hpp"

int main()
{
Atomic_Structure gas;
gas.read_xyz("gas.xyz");
for (i=0;i<gas.Nat;i++)
{
for(k=0;k<3;k++)
{
gas.atom[i].v[k]=random_number(-10,10);
}
}
Simulation_Cell cajita;
cajita.periodicity=true;
gas.molecular_dynamic(cajita,1000,"dinamica_gas.xyz");
//gas.geometry_optimization("opt.xyz",3000);
return 0;
}
