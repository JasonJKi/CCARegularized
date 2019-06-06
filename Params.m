classdef Params < handle
    properties
        d 
        kx
        ky
    end
    
    methods (Access = public)

        function this = Params(kx, ky)
            if nargin < 1
                return
            end
            set(this, kx, ky);            
        end
        function this = this(this)
        end
        function this = set(this, kx, ky)
            this.kx = kx;
            this.ky = ky;
            setMinDim(this, kx, ky)         
        end

        function [kx, ky, d] = get(this)
            kx = this.kx;
            ky = this.ky;
            d = this.d;
        end
        
        function setMinDim(this, kx, ky)
            this.d = min(kx, ky);
        end
        
    end
end