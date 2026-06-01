classdef UtilityMethod<handle
    %% The class for string data for dependency (interactions among variable) of the problem"""
    methods (Static)
        function R=gen_rot_mat_pseudo(uu,vv,alp) % Generates a random rotation matrix
            % Created a random rotation matrix
            % see [1] "On the Rigid Rotation Conept in n-Dimensional Spaces"  by
            % Daniele Mortari (2001) for mathematical formulation
            % given two random vectors uu and vv
            u0=uu';v0=vv';
            D=numel(u0);
            u=u0/norm(u0);
            v=v0-(u'*v0)*u;
            v=v/norm(v);
            R=eye(D)+sin(alp)*(v*u'-u*v')+(cos(alp)-1)*(u*u'+v*v');
        end

        function [Y,minDis]=keep_farthest(X,n)
        % iteratively removes closest solutions from X such that n solutions remain in the end
        % each element of X is in [0,1]
            [N,D]=size(X);
            dis=pdist2(X,X);
            dis2=dis+eye(N)*D^.5*2; % distance to self 
            keepInd=[1];
            candidInd=[2:N];
            for k=2:n
                [tmp,index]=max(min(dis(keepInd,candidInd),[],1));% index of farthest point from candidInd
                keepInd=[keepInd candidInd(index)];
                candidInd(index)=[];
            end % for
            Y=X(keepInd,:);
            minDis=min(min(dis2(keepInd,keepInd)));
        end


        function [Y,countTry]=redist_glob_min(X,Xref,hardNich,disTol)
            % distance to the point Xref shrinks according to the rank of the solution
            % when all solutions are sorted according to their distance to Xref % 
            % disTol: Almost makes sure all solutions are at least disTol far from each other
            % 0 <= hardNich: Controls the nonlinearity of scaling
            maxTry=100; % maixumum number of tries so that each redictributed solution is far from the previously relocated solutions
            [N,D]=size(X);
            countTry=zeros(1,N);
            tauNich=0.2;
            Y=X; % redistributed solutions
            if N>1    
                Vlength=zeros(1,N);
                V=zeros(N,D);
                %Vunit=zeros(N,D);
                for k=1:N
                    V(k,:)=X(k,:)-Xref;
                    Vlength(k)=norm(V(k,:));
                    %Vunit(k,:)=V(k,:)/Vlength(k);
                end
                [~,ind]=sort(Vlength,'ascend');
                rnk(ind)=1:N; % rank of each solution when sorted according to distance to Xref
                rnk=(rnk)/(N).*(1-tauNich)+tauNich;
                targetCoef=rnk.^hardNich; % reduce distance to Xref proportionally to the rank
                for k=1:numel(ind)
                    for tryNo=0:maxTry
                        term1=(maxTry-tryNo)/maxTry;
                        finCoef=targetCoef(ind(k))^term1; % final coefficient, ideally close to targetCoef unless it violates the disTol criterion
                        Y(ind(k),:)= Xref+ finCoef*V(ind(k),:);
                        dis1=pdist2(Y(ind(k),:),Y(ind([1:k-1  k+1:N]),:));
                        countTry(ind(k))=tryNo;
                        if  k==1 || k==N ||  (min(dis1)>=disTol)  
                            break % accept this new location
                        end
                    end
                end
            else
                Y=X;
            end
        end


        function [rpr,pr]=calc_robust_peak_ratio(solutions,values,minima,minimum_val,ftol0)
            % Calculate the robust peak ratio given reported solutions X, their values  f, global minima (x_opt),
            % global minimum value (f_opt), and the range of desired tolerance on the values  
            ftol=[min(ftol0), max(ftol0)];
            n_minima=size(minima,1);
            pr=zeros(1,n_minima);
          
            % find corresponding global minimum for each solution
            dis2opt=pdist2(solutions,minima);
            [~,cor_opt_ind]=min(dis2opt,[],2);
            
            % calculate credit for approximating each global minimum
            for k = 1:n_minima
                ind=find(cor_opt_ind==k);
                if numel(ind)>0
                    best_val=min(values(ind));
                    if ftol(2)==ftol(1)
                        pr(k)= (best_val<=(minimum_val+ftol(1)))*1.0;
                    else
                        term1=log(ftol(2))-log(best_val-minimum_val);
                        term2=log(ftol(2))-log(ftol(1));
                        pr(k)=  max(0,min(1,term1/term2));
                    end
                end
            end
            rpr=mean(pr);
        end % function



    end % methods
end % class
