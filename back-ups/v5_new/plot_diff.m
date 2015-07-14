function plot_diff (type, call, list)
z = [] ; 
x = [] ; 
y = [] ; 

for i = 1:numel(type)
    if strcmp(type(i), 'not present ')
        z = [z; list(i)];
    else
        y = [y; list(i)];
    end
    if strcmp(type(i), 'not present ')  && (strcmp(call(i), 'somatic') || strcmp(call(i), 'germline') ) 
        x = [x; list(i)];
    end
end

min_x = 0 
max_x = max([x;y;z])


subplot(3,1,1)
hist(x,100)
xlim([min_x, max_x])
subplot(3,1,2)
hist(y,100)
xlim([min_x, max_x])


subplot(3,1,3)
hist(z,100)
xlim([min_x, max_x])


