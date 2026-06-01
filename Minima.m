classdef Minima<handle
    %% The class for string data for dependency (interactions among variable) of the problem"""
    properties (SetAccess=?ProblemMM)
        X; % (matrix) global minima of the static problem (or at time step #0 if the problem is dynamic)
        f; % the global minimum value (scalar) 
        hard_GO; % hardness of finding each global minimum from GO perspective        
        range_coeff; % global minima are inside this fraction of each dimensionality, excluding close to bounds regions        
        crowd_basin_ind; % index of the global minimum that other solutions are redistributed wrt
        niche_rad; % N
    end
    methods 
        function output=Minima() % requires functions ID and problem dimensionality 
            output.X=nan; % (matrix) global minima of the static problem (or at time step #0 if the problem is dynamic)
            output.f=nan; % the global minimum value (scalar) 
            output.hard_GO=nan; % hardness of finding each global minimum from GO perspective        
            output.range_coeff=0.9; % global minima are inside this fraction of each dimensionality, excluding close to bounds regions        
            output.crowd_basin_ind=nan; % index of the global minimum that other solutions are redistributed wrt
            output.niche_rad=nan; % Niching radius
        end % function
    end % methods
end % class
