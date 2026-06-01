% Basic functions from which the composite functions can be generated
classdef BasicFun  
    % h_GO controls the difficulty of performing global optimization
    % see https://infinity77.net/global_optimization/test_functions.html for function definition
    methods(Static)   
        function f=evaluate(x,hGO,funID)
            if funID==1 
                f=BasicFun.elliptic(x,hGO);  % lambda=0.5
            elseif funID==2
                f=BasicFun.diffPow(x,hGO);   % lambda=0.5,  
            elseif funID==3
                f=BasicFun.schwefelN02Skew(x,hGO); % lambda=0.5
            elseif funID==4
                f=BasicFun.rosenbrock(x,hGO); % lambda=0.5
            elseif funID==5   
                f=BasicFun.ackleySkew(x,hGO); % lambda=0.2,  
            elseif funID==6  
                f=BasicFun.rastrigin(x,hGO); % lambda=1
            elseif funID==7  
                f=BasicFun.weierstrass(x,hGO);   % lambda=5            
            elseif funID==8 
                f=BasicFun.schwefelN26(x,hGO); % lambda=0.005
            else
                throw('This function is not defined')
            end
        end
        function f=ackleySkew(x,hGO)
            skew=5^hGO;            
            coef=skew*(x>0)+1/skew*(x<=0);
            y=coef.*x;
            term1=sqrt(mean(y.^2));
            term2=mean(cos(2*pi*y));
            f=-20*exp(-0.2*term1)-exp(term2)+20+exp(1);
        end
        function  f=diffPow(x,hGO)
            D=numel(x);
            H=hGO*4;
            p=2+H*(0:D-1)/(D-1);
            f=sum(abs(x).^p).^.5;
        end
        function f=elliptic(x,hGO)
            D=numel(x);
            pow=(0:(D-1))/(D-1)*hGO*3;
            f=10000^(0.5-hGO)*sum((10.^pow.*(x)).^2);
        end
        function f=rosenbrock(x,hGO)
            fsphere=20*sum(x.^2);
            y=x+1;
            term1=100*sum((y(2:end)-y(1:end-1).^2).^2);
            term2=sum((y(1:end-1)-1).^2);
            f0=term1+term2;
            f=hGO*f0+(1-hGO)*fsphere;
        end
        function f=rastrigin(x,hGO)
            base=5;
            A=10*(base^hGO-1)/(base-1);
            f= (x.^2)+(A.*(1-cos(2*pi*x)));
            f= sum(f);
        end
        function f=schwefelN02Skew(x,hGO)
            skew=5^hGO;
            D=numel(x);
            coef=skew*(x>0)+1/skew*(x<=0);
            y=x.*coef;
            f0=0;
            for k=1:D
                f0=f0+ sum(y(1:k))^2;
            end
            f=f0;
        end
        function f=schwefelN26(y,hGO)
            base=5;
            H=(base^hGO-1)/(base-1);
            xstar=420.96874635998202731184436501869; % from matlab vpasolve
            fshift=418.9828872724337062747864351956;
            x=y+xstar;
            p1=sum((-300-x).*(x<-500));
            p2= sum((x>500).*(x-420));
            P=p1+p2;
            g=1.5*P+sum(-x.*sin(abs(x).^.5))+fshift*numel(x);
            f=H*g+1*(1-H)*sum(abs(x-xstar));
        end
        function f=weierstrass(y,hGO)  % f14: Weierstrass function (Modified)    
            base=5;
            H=(base^hGO-1)/(base-1);
            D=numel(y);
            a=0.5; % controls the global basin (a higher a makes problem harder)
            b=3; % a higher value makes it makes it more rugged-default is 3
            k=(0:20)'; % level of optima
            term1= (2*pi*b.^k) * (y+0.5);
            h1=(a.^k)'*cos(term1) ;
            h2=sum(a.^k.*cos(pi*b.^k));
            P= sum(  (abs(y)-.5).^1 .*(abs(y)>.5) );
            f0=sum(h1)-D*h2+P;
            f=(f0*H+(1-H)*sum(abs(y)));
        end
    end
end

        
            