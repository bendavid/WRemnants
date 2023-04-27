from copy import deepcopy

from utilities import logging

logger = logging.child_logger(__name__)

class Datagroup(object):

    def __init__(self, name, members=[], scale=None, selectOp=None, selectOpArgs={}, memberOp=None, rebinOp=None, label=None, color=None):
        self.name = name
        self.scale = scale
        self.label = label
        self.color = color

        self.members = members              # list of narf datasets

        self.selectOp = selectOp            # operation that is applied on all members of the group
        self.selectOpArgs = selectOpArgs    # argments to the selectOp
        self.rebinOp = rebinOp              # operation that is applied on all members of the group for changing histogram bins

        self.memberOp = memberOp            # list of operations that is applied on single members

        self.hists = {}                     # list of histograms processed from narf datasets


    def copy(self, new_name, member_filter=None):
        logger.debug(f"Make a copy of group {self.name}, named {new_name}!")

        x = deepcopy(self)
        x.name = new_name

        if member_filter:
            # Invert the member filter and exclude those members
            x.deleteMembers([m for m in filter(lambda x,f=member_filter: not f(x), x.members)])

        return x

    def addMembers(self, members, member_operations=None):
        if callable(member_operations) or member_operations is None:
            for m in members:
                self.addMember(m, member_operations)
        elif len(member_operations) == len(members):
            for m, o in zip(members, member_operations):
                self.addMember(m, o)
        else:
            raise RuntimeError("'member_operations' has to be a string or a list with the same length as 'members'!")            

    def addMember(self, member, member_operation=None):
        # adds a member to the existing members of a given group
        
        # add member operation
        if self.memberOp is None:
            self.memberOp = [None]*len(self.members)
        
        self.memberOp.append(deepcopy(member_operation))
        self.members.append(member)

    def deleteMembers(self, members):
        for m in members:
            self.deleteMember(m)

    def deleteMember(self, member):
        # deletes a process from the list of members of a given group

        if member not in [m for m in self.members]:
            logger.warning(f"The member {member.name} can not be found in the group {self.name}! Do nothing here.")
            return

        logger.debug(f"Delete member {member.name}!")

        mask = [m != member for m in self.members]
        if self.memberOp is not None:
            self.memberOp = [v for (v, i) in zip(self.memberOp, mask) if i]

        self.members = [m for (m, i) in zip(self.members, mask) if i]
