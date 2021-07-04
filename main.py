import torch as t
from torch.nn import Transformer
from torch.nn import MSELoss
from torch.optim import Adam
print("Start program")

a1 = t.rand((10, 10, 100), requires_grad=True)
a2 = t.rand((10, 10, 100), requires_grad=True)
a3 = t.rand((10, 10, 100), requires_grad=True)

model = Transformer(nhead=10, d_model=100)

optimizer = Adam(model.parameters(), lr=0.001)

out = model(a1, a2)
loss = MSELoss()
loss_val = loss(out, a3)

optimizer.zero_grad()
loss_val.backward()
optimizer.step()


