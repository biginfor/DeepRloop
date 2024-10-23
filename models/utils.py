import torch
import torch.optim as optim


def get_optim(
        model: torch.nn.Module, 
        name: str, 
        **kwargs):
    if name.lower() == 'sgd':
        return optim.SGD(
            model.parameters(), 
            lr = kwargs.get('lr', 1e-2),
            momentum=kwargs.get('momentum', 0.9),
            weight_decay=kwargs.get('weight_decay', 1e-4),
            nesterov=kwargs.get('nesterov', True)
        )
    elif name.lower() == 'adamw':
        return optim.AdamW(
            model.parameters(),
            lr = kwargs.get('lr', 1e-2),
            weight_decay=kwargs.get('weight_decay', 1e-4),
        )
    elif name.lower() == 'adam':
        return optim.Adam(
            model.parameters(),
            lr = kwargs.get('lr', 1e-2),
            weight_decay=kwargs.get('weight_decay', 1e-4),
        )
    elif name.lower() == 'rmsprop':
        return optim.RMSprop(
            model.parameters(), 
            lr = kwargs.get('lr', 1e-2),
            momentum=kwargs.get('momentum', 0.9),
            weight_decay=kwargs.get('weight_decay', 1e-4),
        )
    

def get_loss():
    return torch.nn.BCELoss()