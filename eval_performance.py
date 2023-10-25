import json
import argparse
import numpy as np


def load_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)

    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate elicitation performance given some json file')
    parser.add_argument('--filename', type=str, required=True)

    args = parser.parse_args()
    data = load_data(args.filename)
    reg = np.sum([v[0] for k, v in data.items()])/len(data)
    mcts = np.sum([v[1] for k, v in data.items()])/len(data)
    print('Average mcts performance: %s vs average regular elicitation performance: %s' % (mcts, reg))
    print('MCTS is effective by', 1 - mcts/reg)


