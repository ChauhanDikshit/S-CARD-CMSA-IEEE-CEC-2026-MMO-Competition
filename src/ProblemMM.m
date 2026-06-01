% ********* Class of New Composite Multimodal Optimization Problems for CEC 2026 competition **************
% ********* see example1.py in the upper folder for instrcutions on using this class **************
% code developed by Ali Ahrari
% last update on 10-Feb-2026 by A. A.
% random seed numbers have been generated using MATLAB:
    % seed 0: for 10000 uniform numbers with precision of 10
    % seed 1: for normal numbers with precision of 10
    % seed 2: for 248 permutations of size 10000
    % seed 3: for 8 global minimum values between (-500,500)
 


classdef ProblemMM<handle
    properties (Access=public)
        max_instance_no; % a total of 15 problem instances 
        pid; % problem ID
        fun_id; % function ID for the basic function (different from problem ID) 
        n_minima; % (scalar) number of global minima
        hard_GO; % a number in [0,1] that specifies the hardness from global optimization perspective
        hard_NU; % a non-negative Real number specifying the non-uniformity in the distribution of global minima 
        instance_no; % problem instance No
        dim; % dimensionality
        low_bound; %the lower bound of the search space
        up_bound; % the upper bound of the search space       
        lambda0; % for scaling the search range of thmax_evale basic function
        d_min; % distance threshold between global minima
        max_eval_coeff; % the coefficient for the evaluation budget
        init_mult; % candidate multiplier for samling uniform solutions for coordinates of global minima
        sigma_width;  % controls the impact extent of a basic function  
    end
    properties (SetAccess =protected)
        used_eval; % used evaluation so far   
        master_seed; % problem special number used for reading random numbers from CSV files
        max_eval; % the evaluation budget
        depend; % the dependency structure among variables
        rotation; % data for rotation of the modes
        minima; % information about the global minima (locations,...): For postprocessing of results only
        normal_numbers; % a series of random numbers from the standard normal distribution
        uniform_numbers; % a series of random numbers from the standard uniform distribution
        index_normal; % current index for reading number from normal_numbers  
        index_uniform; % current index for reading number from uniform_numbers  
    end
    methods      
        function problem=ProblemMM(pid,instance_no,dim) % requires functions ID and problem dimensionality         
            problem.max_instance_no=15; % a total of 15 problem instances 
            if (instance_no>problem.max_instance_no) || (instance_no<1) || (floor(instance_no)~=instance_no) 
                error('Error: instance_no should be an integer between 1 and 15 (inclusive)')
            end
            if pid<1 || pid>16 || floor(pid)~=pid
                error('Error: pid should be an integer between 1 and 16 (inclusive)')
            end
            if dim<2 || dim>80 || floor(dim)~=dim
                error('Error: dim should be an integer between 2 and 50 (inclusive)')
            end
            data=readtable("data/pid-data.xlsx");        
            problem.pid=pid; %problem ID
            problem.fun_id=floor(.5+data.fun_id(pid)); % function ID for the basic function (different from problem ID) 
            problem.n_minima=floor(.5+data.n_minima(pid)); % (scalar) number of global minima
            problem.hard_GO=[data.hard_GO_min(pid),data.hard_GO_max(pid)]; % a number in [0,1] that specifies the hardness from global optimization perspective
            problem.hard_NU=data.hard_NU(pid); % a non-negative Real number specifying the non-uniformity in the distribution of global minima 
            problem.instance_no=instance_no; % problem instance No
            problem.dim=dim; % dimensionality
            problem.low_bound=-5*ones(1,problem.dim); %the lower bound of the search space
            problem.up_bound=5*ones(1,problem.dim); % the upper bound of the search space       
            problem.lambda0=data.lambda0(pid); % for scaling the search range of the basic function
            problem.d_min=0.3*problem.dim^.5; % distance threshold between global minima
            problem.max_eval_coeff=data.max_eval_coeff(pid); % the coefficient for the evaluation budget
            problem.init_mult=5; % candidate multiplier for samling uniform solutions for coordinates of global minima
            problem.sigma_width=.5;  % controls the impact extent of a basic function  
    
            % protected variables
            problem.used_eval=0; % used evaluation so far   
            problem.master_seed=nan; % problem special number used for reading random numbers from CSV files
            problem.max_eval=nan; % the evaluation budget
            problem.depend=Dependency(problem.pid,problem.dim,problem.instance_no,problem.max_instance_no); % the dependency structure among variables
            problem.rotation=Rotation(); % data for rotation of the modes
            problem.minima=Minima(); % information about the global minima (locations,...): For postprocessing of results only
            problem.normal_numbers=nan; % a series of random numbers from the standard normal distribution
            problem.uniform_numbers=nan; % a series of random numbers from the standard uniform distribution
            problem.index_normal=1; % current index for reading number from problem.normal_numbers  
            problem.index_uniform=1; % current index for reading number from problem.uniform_numbers  
        end
 
 
        function form(problem) % ******************* formulate problem *******************************
            % ******************* Load problem data from CSV files *********************      
            num_uniform = readmatrix("data/num-uniform.csv"); % array of uniformly distributed random numbers in (0,1)
            num_normal = readmatrix("data/num-normal.csv"); % array of numbers with standard normal distribution 
            sequences = floor( .5 + readmatrix("data/sequences.csv")); % matrix, permutations of the aforementioned random numbers to be used successively
            fstar_data=readmatrix("data/fstar-data.csv"); % global minimum values - 1-D array
    
            % ****************************** set the evaluation budget and the global minimum value **************************** 
            problem.max_eval=floor(.5+problem.max_eval_coeff*problem.dim);
            problem.minima.f=fstar_data(problem.fun_id);

            % ***** specify array of random numbers (uniform and normal) to be used for benchmark generation ******  
            problem.master_seed=problem.max_instance_no*(problem.pid-1)+problem.instance_no; % index number for this pid and instance_no
            used_seq=sequences(problem.master_seed,:); % use this sequence of random numbers
            problem.normal_numbers=num_normal(used_seq); % rearranged random numbers from normal distribution
            problem.uniform_numbers=num_uniform(used_seq); %rearranged random numbers from uniform distribution
            clear used_seq

            problem.det_minima_coords() % specify the locations of global minima
            problem.det_minima_hardness() % determine the hardness of each mode from global optimization or convergence perspective
            problem.det_niche_rad() % the niche radius for each global minimum based on the half of the distance to the closest global minima
            problem.det_depend_base_struct() % form the base dependency structure 
            problem.det_depend_struct() % form the dependency structure for each mode by perturbation of the base dependency struture 
            problem.det_rotation_mat() % create the rotation matrices given the dependency structures


        end % function
        function det_minima_coords(problem)
            % ********************** determine locations of global minima ***********************        
            temp=problem.uniform_numbers(problem.index_uniform:problem.index_uniform+problem.init_mult*problem.dim*problem.n_minima-1); % select random uniform numbers several times of n_minima*dim
            problem.index_uniform=problem.index_uniform+numel(temp); % for reading subsequent numbers 
            rand_points=reshape(temp,problem.dim,problem.n_minima*problem.init_mult)'; % solutions with random distribution from which global minima are selected
            % set the reference solution for redistribution 
            problem.minima.crowd_basin_ind=floor(ceil(problem.uniform_numbers(problem.index_uniform)*problem.n_minima)); % index of Xref is selected randomly
            problem.index_uniform=problem.index_uniform+1; % update the index of used random number with uniform distribution
            % Now create a relatively uniformly distributed set of points from randomly distributed points
            [uniform_points,~]=UtilityMethod.keep_farthest(rand_points,problem.n_minima); % select farthest ones (remove closest ones iteratively)
            uniform_points=uniform_points(1:problem.n_minima,:)*problem.minima.range_coeff+(1-problem.minima.range_coeff)/2; % make sure minima are not too close to the bounds
            uniform_points = uniform_points.*repmat(problem.up_bound-problem.low_bound,problem.n_minima,1)+repmat(problem.low_bound,problem.n_minima,1); % map uniform distribution from [-1,1]^dim to search space
            % Now set global minima locations by redistributing the uniform points (to make them non-uniform)
            [problem.minima.X,~]=UtilityMethod.redist_glob_min(uniform_points,uniform_points(problem.minima.crowd_basin_ind,:),problem.hard_NU,problem.d_min); % non-uniform distribution
        end
        function det_minima_hardness(problem)  
            % ******************************** determine the hardness of each global minimum  *************************** %             
            [~,ind0]=sort(problem.uniform_numbers(problem.index_uniform:problem.index_uniform+problem.n_minima-1)); % use random numbers to sort out the hardness
            problem.index_uniform=problem.index_uniform+problem.n_minima;
            if problem.n_minima>1
                coef=(ind0-1)/(problem.n_minima-1); % This is n_minima uniformly distributed numbers in [0,1] with random order
                problem.minima.hard_GO= coef* (problem.hard_GO(2)-problem.hard_GO(1))+problem.hard_GO(1); % hard_GO for modes uniformly changes from the problem.hard_GO[0] to problem.hard_GO[1]
            else % in an unwanted case when there is only one global mode
                problem.minima.hard_GO=.5*problem.hard_GO(1)+.5*problem.hard_GO(2);
            end
        end
        
        function det_niche_rad(problem)
            % ************************************** determine the niche_rad for each global minimum ******************************* %
            if problem.n_minima==1
                problem.minima.niche_rad=5*sqrt(problem.dim);
            else % there are more than one global minima
                tmp=pdist2(problem.minima.X,problem.minima.X);
                tmp=tmp+max(max(tmp))*eye(problem.n_minima); % ignore diagonal elements (distance to problem)
                problem.minima.niche_rad=min(tmp,[],1)/2.0; % half of distance to the closest global minimum
            end
        end

        function det_depend_base_struct(problem)        
            % ****************************** determine base dependency structure for variables ******************************** %        
            if problem.depend.n_blocks==problem.dim % fully separable
                problem.depend.base_struct=mat2cell((1:problem.dim)', ones(1,problem.dim) ,1 )';
            elseif problem.depend.n_blocks==1 % fully rotated
                problem.depend.base_struct=mat2cell(1:problem.dim,1,problem.dim);
            elseif problem.depend.n_blocks<problem.dim && problem.depend.n_blocks>1 % block separability
                temp=problem.uniform_numbers(problem.index_uniform:problem.index_uniform+problem.dim-1); problem.index_uniform=problem.index_uniform+problem.dim;
                [~,rand_perm]=sort(temp); % a random permutation of dimensions (0 to dim-1)
                problem.depend.block_sizes=problem.depend.block_size_min*ones(1,problem.depend.n_blocks); % minimum size of each block is 1
                candid_blocks=1:problem.depend.n_blocks; % all blocks can increase their size

                while sum(problem.depend.block_sizes)<problem.dim % until the sum of the sizes of all blocks equald to problem dimensionality
                    temp=problem.uniform_numbers(problem.index_uniform); problem.index_uniform=problem.index_uniform+1;
                    ind=ceil(temp*numel(candid_blocks)); % choose of the candidate blocks randomly
                    ind2=candid_blocks(ind); % index of block to increase in size
                    problem.depend.block_sizes(ind2)=problem.depend.block_sizes(ind2)+1; % increase the size of this block by one
                    if problem.depend.block_sizes(ind2)>=problem.depend.block_size_max % if this block size is equal or greater than the upper limit
                        candid_blocks=setdiff(candid_blocks, ind2); % do not consider this block for enlarging in the next iteration
                    end
                end % while
                % now given the block sizes and random ordering of variables (rand_perm), assigns variables to blocks 
                ind=1;
                for k = 1:problem.depend.n_blocks
                    these_dims=rand_perm(ind:ind+problem.depend.block_sizes(k)-1);
                    problem.depend.base_struct{k}=these_dims;
                    ind=ind+problem.depend.block_sizes(k);
                end
            end % if
        end

        function det_depend_struct(problem)         
            % ************* set all dependency structures for all modes by perturbation of base dependency structure ************** %        
            problem.depend.struct=cell(1,problem.n_minima);
            for glob_no = 1:problem.n_minima
                problem.depend.struct{glob_no} =  problem.depend.base_struct;  % the dependency structure of the mode initially gets the base dependency structure 
                if problem.depend.n_blocks<problem.dim && problem.depend.n_blocks>1 % apply random perturbation (unless special case of fully separable or fully rotated, for each perturbation is meaningless)
                    for swap_count = 1:problem.depend.n_swap % apply a predefined number of swaps 
                        temp=problem.uniform_numbers(problem.index_uniform:problem.index_uniform+2-1); problem.index_uniform=problem.index_uniform+numel(temp);
                        indexes= ceil(temp*problem.dim);% choose two indexes from 0 to dim-1 to be swapped
                        % find the block_no and element number for each
                        block_size_cumsum=cumsum(problem.depend.block_sizes); % cumulative sizes of block sizes
                        block_ind = 1+ sum(reshape(indexes,[],1)>block_size_cumsum,2); % indexes of blocks of dimensions to swap
                        block_size_cumsum_plus=[0 block_size_cumsum]; % put zero before cumulative sum
                        index_in_block=indexes-block_size_cumsum_plus(block_ind); % indexes in the corresponding blocks
                        % Now having indexes of blocks and indexes at blocks of both variables, swap them 
                        temp=problem.depend.struct{glob_no}{block_ind(1)}(index_in_block(1))+0; % keep the first value that should be swapped
                        problem.depend.struct{glob_no}{block_ind(1)}(index_in_block(1))=problem.depend.struct{glob_no}{block_ind(2)}(index_in_block(2))+0;
                        problem.depend.struct{glob_no}{block_ind(2)}(index_in_block(2))=temp;
                    end
                end
                if 1 % c=heck point only-
                    term1=[];
                    term2=[];
                    for k = 1:problem.depend.n_blocks
                        term1=[term1 problem.depend.base_struct{k}];
                        term2=[term2,problem.depend.struct{glob_no}{k} ];
                    end
                    if sum(term1==term2)<problem.dim-problem.depend.n_swap*2  % if variation is more than the upper limit
                        error('Error: variation in dependency blocks is more than the upper limit')
                    end
                    if max(abs(sort(term1)-(1:problem.dim)))>.0001
                        error('Error: some dimensions are not in the base dependency structure')
                    end
                    if max(abs(sort(term2)-(1:problem.dim)))>.0001
                        error('Error: some dimensions are not in the dependency strcuture for this mode')
                    end
                end % checkpoint if 
            end % for each glob_no
        end % function

        function det_rotation_mat(problem)
            % ******************************* determine the rotation matrices for modes*********************************** %
            
            problem.rotation.mat=cell(1,problem.n_minima); % pre-allocation
            for glob_no = 1:problem.n_minima % for each mode
                problem.rotation.mat{glob_no}=eye(problem.dim);
                keep=1:problem.dim; % for checkpoint
                for block_no = 1:problem.depend.n_blocks % for each block of this mode
                    block_size=numel(problem.depend.struct{glob_no}{block_no}); 
                    if block_size>1 % create the subspace rotation matrix
                        temp0=2*block_size; % you need this number of random numbers (normal distribution) to form the rotation matrix for this block
                        temp_uv=problem.normal_numbers(problem.index_normal:problem.index_normal+temp0-1); problem.index_normal=problem.index_normal+temp0; % get first temp0 numbers of the sequence of Normal numbers
                        u0=temp_uv(1:block_size); % First random vector to create the rotation matrix
                        v0=temp_uv(block_size+1:end); % % Second random vector to create the rotation matrix
                        angle=problem.uniform_numbers(problem.index_uniform)*problem.rotation.angle_max; problem.index_uniform=problem.index_uniform+1; % the rotation angle
                        subspace_rot_mat=UtilityMethod.gen_rot_mat_pseudo(u0,v0,angle); % create the subspace rotation matrix 
                        % now replace the corresponding elements of the full rotation matrix for this mode by the elements of the subspace rotation matrix
                        temp=problem.depend.struct{glob_no}{block_no};                    
                        problem.rotation.mat{glob_no}(temp,temp)=subspace_rot_mat(:,:);  
    
                        if ~all(ismember(temp, keep)) % checkpoint: the modified rows/columns must not have previously been modified because there is no overlap among blocks
                            error('Error: Why there is an overlap between subspaces for rotation?')
                        end
                        keep=setdiff(keep,temp);
                    end % subspace rotation matrix was created and included
                end % for each block
            end % for each mode
        end% function

        % objective function accepts a matrix where each row is a solution
        function f=func_eval(problem,x) 
            N=size(x,1);
            f=zeros(1,N);
            for k = 1:N
                f(k)=problem.func_eval_single(x(k,:))+problem.minima.f;
            end
        end
             
    
        function f=func_eval_single(problem,x) % calculate the fitness of the solution x
            F=zeros(1,problem.n_minima); % fitness values
            for k = 1:problem.n_minima % calculate the fitness value for each basic function independently
                shift=problem.minima.X(k,:);
                x_rot=(x-shift) * problem.rotation.mat{k};
                F(k)=BasicFun.evaluate(x_rot/problem.lambda0,problem.minima.hard_GO(k),problem.fun_id);
            end
            % Calculate weights
            dis=pdist2(x,problem.minima.X); % distance of the point to all global basins
            norm_dis2=(dis./(problem.sigma_width*problem.minima.niche_rad)).^2;
            norm_dis2min=min(norm_dis2);
            if norm_dis2min<=1 % some addition in case all distances are too large and all weights might become zero
                C0=0;
            else
                C0=1-norm_dis2min;
            end
            weights=exp(-norm_dis2-C0); % raw weight of each basic function on the fitness of solution x
            max_weights=max(weights); % highest weight
            chk= abs(weights-max_weights)<1e-14;
            weights=weights.*(1-max_weights.^10) .* (1-chk) + weights .* chk; % the highest raw weight does not change, the rest share the rest
            weights=weights/sum(weights); % adjusts weights 
            f=sum(weights*F');  % the  fitness is the weighted average all basic functions        
            problem.used_eval=problem.used_eval+1; % update this property
        end % function   

 
    end % methods
end % class

 

