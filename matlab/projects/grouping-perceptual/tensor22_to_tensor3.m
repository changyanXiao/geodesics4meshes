function T = tensor22_to_tensor3(U)

T = cat(3, U(:,:,1,1), U(:,:,2,2), U(:,:,1,2));

end