function [newarr] = inrange(arr,range)
    newarr=double(arr >= range(1) & arr <= range(2));
end