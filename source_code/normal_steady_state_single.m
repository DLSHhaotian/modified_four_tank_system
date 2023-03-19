function s_yu=normal_steady_state_single(y,u,ys,us,index_step)
[y_row,y_col]=size(y);
s_yu=zeros(y_row,y_col);
for i=1:y_col
    s_yu(i)=(y(i)-ys)/(u(i)-us);
end
s_yu(1:index_step-1)=s_yu(index_step);
end