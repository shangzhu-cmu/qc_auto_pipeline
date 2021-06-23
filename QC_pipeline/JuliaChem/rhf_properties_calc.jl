#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem
import Statistics
#import HypothesisTests 
import MPI
using BenchmarkTools

function execute(molecule, model, keywords)
  mol, basis = JuliaChem.JCBasis.run(molecule, model; output=2)
  scf = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"]; 
    output=2)
  prop = JuliaChem.JCRHF.Properties.run(mol, basis, scf, keywords["prop"];
    output=2) 
end

#================================#
#== JuliaChem execution script ==#
#================================#
function script(input_file)
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  try
    #== perform scf benchmark ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file; 
        output=0)
 
    if (driver == "energy")
      if (model["method"] == "RHF")
        timings = BenchmarkTools.@benchmark execute($molecule, $model, $keywords) 
        display(timings)
      end
    end
    MPI.Barrier(MPI.COMM_WORLD)
    
    #== finalize JuliaChem runtime ==#
    JuliaChem.finalize()
  
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end    
end

script(ARGS[1])
