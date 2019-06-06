classdef CCA < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        params = Params();
        covMatrix = CovMatrix();
        
        % input data
        x
        y
        
        % canonical component weights
        A
        B
        
        % canonical components
        U
        V
    end
    
    methods
        
        function this = CCA(params)
            this.params = params.this();          
        end
        
        function computeCov(this, x, y)
            [rxx, ryy, rxy, ryx] = this.computeNanCov(x, y);
            this.covMatrix.set(rxx, ryy, rxy, ryx)
        end
        
        function [x, y] = normalize(this, x, y)
            x = this.nanMeanNormalization(x);
            y = this.nanMeanNormalization(y);
        end
        
        function [A, B] = computeMaximizedCorrComponents(this, covMatrix, params)
            [kx, ky, d] = params.get();

            [rxx, ryy, rxy, ryx] = covMatrix.get();
            
            [rxxRegSqrtInvs, ryyRegInvs, ryyRegSqrtInvs] = this.regularizedInverseCov(rxx, ryy, kx, ky);
            
            M =  rxxRegSqrtInvs * rxy * ryyRegInvs * ryx * rxxRegSqrtInvs; 
            M = (M+M') / 2;
            
            [C,Dc] = eig(M);  
            [~,indx ] =sort(diag(Dc), 'descend');
            C = C(:, indx(1:d)); 
            A = rxxRegSqrtInvs * C;
            
            D = ryyRegSqrtInvs * ryx * rxxRegSqrtInvs * C;
            B = ryyRegSqrtInvs * D;
        end
        
        function [rho, p] = predict(this, x, y, A, B)
            if nargin < 4
                A = this.A;
                B = this.B;
            end
            
            x = this.noNan(x);
            y = this.noNan(y);
            [U, V] = this.computeComponents(x, y, A, B);
            
            setComponents(this, U, V);
            
            [rho, p] = computeCorrelation(U, V);
        end
        

        
        function [U, V] = getComponents(this)
            U = this.U;
            V = this.V;
        end
        
        function setComponents(this, U, V)
            this.U = U;
            this.V = V;
        end
        
        
        function this = svdRegularize(this, str)
            [kx, ky] = getHyperParams(this);
            switch str
                case 'x'
                    rxx = this.covMatrix.rxx;
                    kx = this.numComponentsExplainingEigenValVariance(rxx);
                case 'y'
                    ryy = this.covMatrix.ryy;
                    ky = this.numComponentsExplainingEigenValVariance(ryy);
            end
            setHyperParams(this, kx, ky)
        end            
        
        function setWeights(A, B)
            this.A = A;
            this.B = B;
        end
        
        function [A, B] = getWeights(this)
            A = this.A;
            B = this.B;
        end
        
        function fit(this, x, y)
             % Normalize train data and assign zero mean (while ignoring nan values).
            [x, y] = normalize(this, x, y);           
            computeCov(this, x, y)
            [this.A, this.B] = computeMaximizedCorrComponents(this, this.covMatrix, this.params);
        end
        
    end
    
    methods (Static)
        
        function [rxx, ryy, rxy, ryx]  = computeNanCov(x, y)
            [rxx, ryy, rxy, ryx] = nanRXY(x, y);
        end
        
        function [U, V] = computeComponents(x, y, A, B)
            U = x*A;
            V = y*B;
        end
        
        
        function [rxxRegSqrtInvs, ryyRegInvs, ryyRegSqrtInvs] = regularizedInverseCov(rxx, ryy, kx, ky)
            rxxRegSqrtInvs = regSqrtInv(rxx, kx); % regularized Rxx^(-1/2)
            ryyRegInvs = regInv(ryy, ky);  % regularized Ryy^-1
            ryyRegSqrtInvs = regSqrtInv(ryy, ky); % regularized Ryy^(-1/2)
        end
        
        function numComponent = numComponentsExplainingEigenValVariance(rxx, eigValThresh)
            if nargin < 2
                eigValThresh = 0.99;
            end
            [v, d] = eig(rxx);
            numComponent = find(cumsum(diag(d)/sum(diag(d))) > eigValThresh,1);
        end
       
        
        function x = nanMeanNormalization(x)
            n = size(x,1); % n - number of samples
            x = x - repmat(nanmean(x,1), n, 1);
        end
        
        function x = noNan(x)
            nanInd = isnan(x);
            x(nanInd) = 0;
        end
        
    end
end