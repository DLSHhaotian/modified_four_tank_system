function s_yu=normal_steady_state(y,u,ys,us,index_step,index_u)
[y_row,y_col]=size(y);
s_yu=zeros(y_row,y_col);
for i=1:y_col
    for j=1:4
        s_yu(j,i)=(y(j,i)-ys(j))/(u(index_u,i)-us(index_u));
    end
end
s_yu(1,1:index_step-1)=s_yu(1,index_step);
s_yu(2,1:index_step-1)=s_yu(2,index_step);
s_yu(3,1:index_step-1)=s_yu(3,index_step);
s_yu(4,1:index_step-1)=s_yu(4,index_step);
end