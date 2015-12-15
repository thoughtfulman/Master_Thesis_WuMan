classdef ChannelEstiMethod < Simulink.IntEnumType
       enumeration
           % enum channel estimation method
            IDEAL (1)
            LS (2)
            LMMSE (3)
            SVD (4) 
       end
end