import { expect as expectCDK, haveResource } from '@aws-cdk/assert';
import * as cdk from '@aws-cdk/core';
import { core } from '@myhelix/cdk-library';
import * as Cdk from '../lib/cdk-stack';
import { Local } from "../lib/local";

test('stack has components', () => {
    const app = new cdk.App();
    const local = new Local();
    // WHEN
    const stack = new Cdk.R2vStack(app, 'MyTestStack', {
        accountingTag: core.AccountingCategory.ENGINEERING,
        serviceTag: `test-${local.projectName}-cdk`,
        namedEnv: core.Environment.PlatformDevelopment
    });
    // THEN
    expectCDK(stack).to(haveResource('AWS::Batch::JobDefinition'));
});
