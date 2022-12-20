

import torch


class BaseModel(torch.nn.Module):

    def __init__(self, hyper_parameters):
        """

        :param hyper_parameters: {}
        """
        super(BaseModel, self).__init__()

        self.hyper_parameters = hyper_parameters
        inp_s=self.hyper_parameters['inp_shape']
        self.linear1 = torch.nn.Linear(self.hyper_parameters['inp_shape'], inp_s*2)

        self.activation = torch.nn.ReLU()
        self.linear2 = torch.nn.Linear(inp_s*2, inp_s*4)
        self.linear3 = torch.nn.Linear(inp_s*4, inp_s*4)
        self.linear4 = torch.nn.Linear(inp_s * 4, inp_s * 2)
        self.linear5 = torch.nn.Linear(inp_s * 2, self.hyper_parameters['output_shape'])

        # self.softmax = torch.nn.Softmax(dim=1)

    def forward(self, x):
        # print("x_shape:", x.shape)
        x = self.linear1(x)
        x = self.activation(x)
        x = self.linear2(x)
        x = self.activation(x)
        x = self.linear3(x)
        x = self.activation(x)
        x = self.linear4(x)
        x = self.activation(x)
        x = self.linear5(x)
        x = self.activation(x)
        return x
