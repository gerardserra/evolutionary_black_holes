
function weights = stepWeights05(total, step)
    weights = ones(total, 1);
    for i = 1:total
        if i < step
            weights(i) = weights(i) - ((1-0.5)/step)*i;
        else
            weights(i) = 0.5;
        end
    end
end