function [ Tdata, Tflux ] = ...
    pmMorph( Adata, ...      % start photon position 2xN
             Idata, ...   % position after space-deformation
             Bdata, ...      % end photon position 2xN
             permutation, ...   % end points of a match
             Aflux, ...      % start flux
             Bflux, ...      % end flux
             blend, ...      % (1-blend)*A + blend*B
             ratweight)         % weight for the rational quadratic blending
%PMMORPH Morphs point sets into each other based on assigned matches

    
    %% blending the position
    Bdata = Bdata(:,permutation);
    Bflux = Bflux(:,permutation);
	
    if size(Idata,1) == 0 || ratweight == 0
        Tdata = Adata + (Bdata-Adata) .* blend;
    else
        Tdata = pmInterpolatePnts(Adata, Idata, Bdata, blend, ratweight);
    end
    
    %% blend the flux
    Tflux = Aflux + (Bflux-Aflux) .* blend;

end

