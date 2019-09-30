function [val_min, pos_min]=min_array(A)
    val_min=A(1,1);
    for i=1:length(A)
        if A(1,i)<val_min
            val_min=A(1,i);
            pos_min=i;
        end
    end
end
