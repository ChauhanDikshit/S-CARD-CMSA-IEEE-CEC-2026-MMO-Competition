classdef Dependency<handle
    %% The class for string data for dependency (interactions among variable) of the problem"""
    properties (SetAccess=?ProblemMM)         
        n_blocks; % number of blocks
        block_size_delta_coeff; % controls variation in 
        n_swap; % number of random swaps applied to each connectivity 
        block_size_mean; % mean value for the block size
    end
    properties (SetAccess=?ProblemMM)
        block_size_min; % minimum value for block size
        block_size_max;
        block_sizes;  % actual sizes for the blocks      
        base_struct;     % mean (base) dependency structure 
        struct=nan;
    end
    methods 
        function depend=Dependency(pid,dim,instance_no,max_instance_no) % requires functions ID and problem dimensionality 
            data=readtable("data/pid-data.xlsx");        

            depend.n_blocks=1+floor((instance_no-1)/(max_instance_no-1)*(dim-1)+.5);
            depend.block_size_delta_coeff=data.block_size_delta_coeff(pid);
            depend.n_swap=1+floor(dim*data.n_swap_coeff(pid)); % number of random swaps applied to each connectivity 
            depend.block_size_mean=1/depend.n_blocks*dim; % mean value for the block size
            temp=floor(min([floor(depend.block_size_mean) , floor(.5+depend.block_size_mean*(1-depend.block_size_delta_coeff))] ));
            depend.block_size_min=max([1,temp]); % minimum value for block size
            depend.block_size_max=floor(max([ceil(depend.block_size_mean) , floor(.5+depend.block_size_mean*(1+depend.block_size_delta_coeff))]));
            depend.block_sizes=nan;  % actual sizes for the blocks      
            depend.base_struct=cell(1,depend.n_blocks);     % mean (base) dependency structure 
            depend.struct=nan;
        end % function
    end % methods
end % class
