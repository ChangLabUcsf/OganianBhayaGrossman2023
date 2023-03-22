function [cols] = getColors(type)
    switch type
        case 1 % vowels
            cols = brewermap(5, 'Dark2');
            cols(5, :) = [0.2 0.6 0.9];
    end
end