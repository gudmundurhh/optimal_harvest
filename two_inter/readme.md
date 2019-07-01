Computing the optimal harvesting policy for the base model is done in the script hjb_2dMAIN.py. 
It uses the boundary value functions which are described in boundary_value.py.
hjb_2dMAIN.py is also identical to hjb_base.py, the only difference being that the latter is a callable fucntion that is used when comparing the value function and the optimal policies for the base model to their counterparts in the cases.

The optimal harvesting policy for the adjusted noise system is generated in thehjb_noise.py. There it is shown how the value functions are compared. 
