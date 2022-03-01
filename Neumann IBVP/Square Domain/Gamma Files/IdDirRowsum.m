function [F] = IdDirRowsum(EdgeTab,Tj)
    A = zeros(Tj,1);

        F = find(~EdgeTab); %Finds index of the 0 entries in EdgeTab
        if isempty(F) == 1 %if F is a empty vector (i.e. no 0 entries in EdgeTab), we force F to be 0
            F=0;
        end

    F = [F; zeros(numel(A)-numel(F),1)]; %Constructs vector of size Tj, with the index included.  
     
    
end
