function val=data_point_estimate(X,Y,intervals)
%%%%%%%%%The function was used to estimate the value in y axis
n=length(intervals);
m=length(X);

for i=1:n
    for j=1:m
        if (intervals(i) < X(j))
            break;
        end
    end
    
    if (intervals(i)==X(j-1))
        val(i,1)=Y(j-1,1);
    else
        f=fit([X(j-1,1);X(j,1)],[Y(j-1,1);Y(j,1)],'poly1');
        val(i,1) = f.p1*intervals(i) + f.p2;
    end
    
end

end