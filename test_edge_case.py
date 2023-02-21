from rdkit import Chem


mol = Chem.MolFromSmiles("O(c1ccc(cc1)CNc2nc(nc3N(C=Nc32)C(C)C)N)C")
purine_true = Chem.MolFromSmiles("c1c2c(nc[nH]2)ncn1")
purine_fake = Chem.MolFromSmiles("c1c2c(nc[NH]2)ncn1")
old_smarts = Chem.MolFromSmarts("c1nc2c([nX3]1)ncnc2") # X should describe all bonds and this nitrogen should have two to the ring and one outside

# But here, when the nitrogen has one hydrogen it doesnt match
print(mol.HasSubstructMatch(old_smarts))
print(purine_true.HasSubstructMatch(old_smarts))
print(purine_fake.HasSubstructMatch(old_smarts))

smarts = Chem.MolFromSmarts("c1nc2c([nv3]1)ncnc2") # This smarts fixes it which describes bond order rather than number of bonds


# So i looked if this is also true for other ringsystems.... and it isnt...
pyrazole_smarts = Chem.MolFromSmarts("c1[nX2][nX3]cc1")
pyarzole_true = Chem.MolFromSmiles("c1cn[nH]c1")
pyrazole_subst = Chem.MolFromSmiles("c1cnn(CCC)c1")

print(pyarzole_true.HasSubstructMatch(pyrazole_smarts))
print(pyrazole_subst.HasSubstructMatch(pyrazole_smarts))

# What dont I understand right now?!