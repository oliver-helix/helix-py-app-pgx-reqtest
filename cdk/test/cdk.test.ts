import { expect as expectCDK, haveResource } from '@aws-cdk/assert';
import { App } from '@aws-cdk/core';
import { core } from '@myhelix/cdk-library';
import { PGXStack } from '../lib/pgx-stack';
import { Local } from "../lib/local";

test('stack has components', () => {
    const app = new App();
    const local = new Local();
    // WHEN
    const stack = new PGXStack(app, 'MyTestStack', {
        accountingTag: core.AccountingCategory.ENGINEERING,
        serviceTag: `test-${local.projectName}-cdk`,
        namedEnv: core.Environment.PlatformDevelopment
    });
    // THEN
    expectCDK(stack).to(haveResource('AWS::Batch::JobDefinition'));
});
