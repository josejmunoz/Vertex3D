function smoothed_data = weighted_moving_average(data, weights, cutoff)
    smoothed_data = [];
    weights_normalised = weights/sum(weights);
    for numElem = 1:length(data)
        valueElem = 0;
        for weight_elem = 1:length(weights_normalised)
            correctedIndex = numElem - (length(weights) - weight_elem);
            if correctedIndex > 0
                valueElem = valueElem + weights_normalised(weight_elem) * data(correctedIndex);
            end
        end
        smoothed_data(numElem) = valueElem;
    end
    smoothed_data(smoothed_data > cutoff) = cutoff;
end

